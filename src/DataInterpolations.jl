module DataInterpolations

### Interface Functionality

abstract type AbstractInterpolation{T, N} end

using LinearAlgebra, RecipesBase
using PrettyTables
using ForwardDiff
using EnumX
import FindFirstFunctions: searchsortedfirstcorrelated, searchsortedlastcorrelated,
                           Guesser

@enumx ExtrapolationType None Constant Linear Extension Periodic Reflective

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
    @inbounds for i in eachindex(t)
        u[:, i] = interp(t[i])
    end
    u
end
function (interp::AbstractInterpolation)(u::AbstractVector, t::AbstractVector)
    @inbounds for i in eachindex(u, t)
        u[i] = interp(t[i])
    end
    u
end

const EXTRAPOLATION_ERROR = "Cannot extrapolate as `extrapolate` keyword passed was `false`"
struct ExtrapolationError <: Exception end
function Base.showerror(io::IO, ::ExtrapolationError)
    print(io, EXTRAPOLATION_ERROR)
end

const LEFT_EXTRAPOLATION_ERROR = "Cannot extrapolate for t < first(A.t) as the `extrapolation_left` kwarg passed was `ExtrapolationType.None`"
struct LeftExtrapolationError <: Exception end
function Base.showerror(io::IO, ::LeftExtrapolationError)
    print(io, LEFT_EXTRAPOLATION_ERROR)
end

const RIGHT_EXTRAPOLATION_ERROR = "Cannot extrapolate for t > last(A.t) as the `extrapolation_right` kwarg passed was `ExtrapolationType.None`"
struct RightExtrapolationError <: Exception end
function Base.showerror(io::IO, ::RightExtrapolationError)
    print(io, RIGHT_EXTRAPOLATION_ERROR)
end

const INTEGRAL_NOT_FOUND_ERROR = "Cannot integrate it analytically. Please use Numerical Integration methods."
struct IntegralNotFoundError <: Exception end
function Base.showerror(io::IO, ::IntegralNotFoundError)
    print(io, INTEGRAL_NOT_FOUND_ERROR)
end

const DERIVATIVE_NOT_FOUND_ERROR = "Derivatives greater than second order is not supported."
struct DerivativeNotFoundError <: Exception end
function Base.showerror(io::IO, ::DerivativeNotFoundError)
    print(io, DERIVATIVE_NOT_FOUND_ERROR)
end

const INTEGRAL_INVERSE_NOT_FOUND_ERROR = "Cannot invert the integral analytically. Please use Numerical methods."
struct IntegralInverseNotFoundError <: Exception end
function Base.showerror(io::IO, ::IntegralInverseNotFoundError)
    print(io, INTEGRAL_INVERSE_NOT_FOUND_ERROR)
end

const INTEGRAL_NOT_INVERTIBLE_ERROR = "The Interpolation is not positive everywhere so its integral is not invertible."
struct IntegralNotInvertibleError <: Exception end
function Base.showerror(io::IO, ::IntegralNotInvertibleError)
    print(io, INTEGRAL_NOT_INVERTIBLE_ERROR)
end

const EXTRAPOLATION_NOT_IMPLEMENTED_ERROR = "The provided extrapolation option is not implemented."
struct ExtrapolationNotImplementedError <: Exception end
function Base.showerror(io::IO, ::ExtrapolationNotImplementedError)
    print(io, EXTRAPOLATION_NOT_IMPLEMENTED_ERROR)
end

export LinearInterpolation, QuadraticInterpolation, LagrangeInterpolation,
       AkimaInterpolation, ConstantInterpolation, QuadraticSpline, CubicSpline,
       BSplineInterpolation, BSplineApprox, CubicHermiteSpline, PCHIPInterpolation,
       QuinticHermiteSpline, LinearInterpolationIntInv, ConstantInterpolationIntInv,
       ExtrapolationType

# added for RegularizationSmooth, JJS 11/27/21
### Regularization data smoothing and interpolation
struct RegularizationSmooth{uType, tType, T, T2, N, ITP <: AbstractInterpolation{T, N}} <:
       AbstractInterpolation{T, N}
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
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
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
            extrapolation_left,
            extrapolation_right)
        N = get_output_dim(u)
        new{typeof(u), typeof(t), eltype(u), typeof(λ), N, typeof(Aitp)}(
            u,
            û,
            t,
            t̂,
            wls,
            wr,
            d,
            λ,
            alg,
            Aitp,
            extrapolation_left,
            extrapolation_right)
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
    T,
    N
} <: AbstractInterpolation{T, N}
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
        N = get_output_dim(u)
        new{typeof(u), typeof(t), typeof(m),
            typeof(p0), typeof(ub), typeof(lb),
            typeof(alg), typeof(pmin), eltype(u), N}(u,
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
