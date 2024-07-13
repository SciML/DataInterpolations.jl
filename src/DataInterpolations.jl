module DataInterpolations

### Interface Functionality

abstract type AbstractInterpolation{T} end

using LinearAlgebra, RecipesBase
using PrettyTables
using ForwardDiff
using ReadOnlyArrays
import FindFirstFunctions: searchsortedfirstcorrelated, searchsortedlastcorrelated,
                           bracketstrictlymontonic

include("parameter_caches.jl")
include("interpolation_caches.jl")
include("interpolation_utils.jl")
include("interpolation_methods.jl")
include("plot_rec.jl")
include("derivatives.jl")
include("integrals.jl")
include("integral_inverses.jl")
include("online.jl")
include("show.jl")

(interp::AbstractInterpolation)(t::Number) = _interpolate(interp, t)
(interp::AbstractInterpolation)(t::Number, i::Integer) = _interpolate(interp, t, i)
function (interp::AbstractInterpolation)(t::AbstractVector)
    u = get_u(interp.u, t)
    interp(u, t)
end

function get_u(u::AbstractVector, t)
    return similar(t, promote_type(eltype(u), eltype(t)))
end

function get_u(u::AbstractVector{<:AbstractVector}, t)
    type = promote_type(eltype(eltype(u)), eltype(t))
    return [zeros(type, length(first(u))) for _ in eachindex(t)]
end

function get_u(u::AbstractMatrix, t)
    type = promote_type(eltype(u), eltype(t))
    return zeros(type, (size(u, 1), length(t)))
end

function (interp::AbstractInterpolation)(u::AbstractMatrix, t::AbstractVector)
    iguess = firstindex(interp.t)
    @inbounds for i in eachindex(t)
        u[:, i], iguess = interp(t[i], iguess)
    end
    u
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

const INTEGRAL_INVERSE_NOT_FOUND_ERROR = "Cannot invert the integral analytically. Please use Numerical methods."
struct IntegralInverseNotFoundError <: Exception end
function Base.showerror(io::IO, e::IntegralInverseNotFoundError)
    print(io, INTEGRAL_INVERSE_NOT_FOUND_ERROR)
end

const INTEGRAL_NOT_INVERTIBLE_ERROR = "The Interpolation is not positive everywhere so its integral is not invertible."
struct IntegralNotInvertibleError <: Exception end
function Base.showerror(io::IO, e::IntegralNotInvertibleError)
    print(io, INTEGRAL_NOT_INVERTIBLE_ERROR)
end

const MUST_COPY_ERROR = "A copy must be made of u, t to filter missing data"
struct MustCopyError <: Exception end
function Base.showerror(io::IO, e::MustCopyError)
    print(io, MUST_COPY_ERROR)
end

export LinearInterpolation, QuadraticInterpolation, LagrangeInterpolation,
       AkimaInterpolation, ConstantInterpolation, QuadraticSpline, CubicSpline,
       BSplineInterpolation, BSplineApprox, CubicHermiteSpline,
       QuinticHermiteSpline, LinearInterpolationIntInv, ConstantInterpolationIntInv

# added for RegularizationSmooth, JJS 11/27/21
### Regularization data smoothing and interpolation
struct RegularizationSmooth{uType, tType, T, T2, ITP <: AbstractInterpolation{T}} <:
       AbstractInterpolation{T}
    u::uType
    û::uType
    t::tType
    t̂::tType
    wls::uType
    wr::uType
    d::Int       # derivative degree used to calculate the roughness
    λ::T2        # regularization parameter
    alg::Symbol  # how to determine λ: `:fixed`, `:gcv_svd`, `:gcv_tr`, `L_curve`
    Aitp::ITP
    extrapolate::Bool
    function RegularizationSmooth(u,
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
        new{typeof(u), typeof(t), eltype(u), typeof(λ), typeof(Aitp)}(
            readonly_wrap(u),
            readonly_wrap(û),
            readonly_wrap(t),
            readonly_wrap(t̂),
            readonly_wrap(oftype(u.parent, wls)),
            readonly_wrap(oftype(u.parent, wr)),
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
    T
} <: AbstractInterpolation{T}
    u::uType
    t::tType
    m::mType        # model type
    p0::p0Type      # initial params
    ub::ubType      # upper bound of params
    lb::lbType      # lower bound of params
    alg::algType    # alg to optimize cost function
    pmin::pminType  # optimized params
    extrapolate::Bool
    function CurvefitCache(u, t, m, p0, ub, lb, alg, pmin, extrapolate)
        new{typeof(u), typeof(t), typeof(m),
            typeof(p0), typeof(ub), typeof(lb),
            typeof(alg), typeof(pmin), eltype(u)}(u,
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

end # module
