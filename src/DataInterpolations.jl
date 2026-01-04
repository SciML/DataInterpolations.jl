module DataInterpolations

### Interface Functionality

abstract type AbstractInterpolation{T} end

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
    if interp.u isa AbstractVector
        # Return a vector of interpolated values, on for each element in `t`
        return map(interp, t)
    elseif interp.u isa AbstractArray
        # Stack interpolated values if `u` was stored in matrix/... form
        return stack(interp, t)
    end
end

function (interp::AbstractInterpolation)(out::AbstractVector, t::AbstractVector)
    if length(out) != length(t)
        throw(DimensionMismatch("number of evaluation points and length of the result vector must be equal"))
    end
    map!(interp, out, t)
    return out
end
function (interp::AbstractInterpolation)(out::AbstractArray, t::AbstractVector)
    if size(out, ndims(out)) != length(t)
        throw(DimensionMismatch("number of evaluation points and last dimension of the result array must be equal"))
    end
    map!(interp, eachslice(out; dims = ndims(out)), t)
    return out
end

const EXTRAPOLATION_ERROR = "Cannot extrapolate as `extrapolate` keyword passed was `false`"
struct ExtrapolationError <: Exception end
function Base.showerror(io::IO, ::ExtrapolationError)
    return print(io, EXTRAPOLATION_ERROR)
end

const LEFT_EXTRAPOLATION_ERROR = "Cannot extrapolate for t < first(A.t) as the `extrapolation_left` kwarg passed was `ExtrapolationType.None`"
struct LeftExtrapolationError <: Exception end
function Base.showerror(io::IO, ::LeftExtrapolationError)
    return print(io, LEFT_EXTRAPOLATION_ERROR)
end

const RIGHT_EXTRAPOLATION_ERROR = "Cannot extrapolate for t > last(A.t) as the `extrapolation_right` kwarg passed was `ExtrapolationType.None`"
struct RightExtrapolationError <: Exception end
function Base.showerror(io::IO, ::RightExtrapolationError)
    return print(io, RIGHT_EXTRAPOLATION_ERROR)
end

const INTEGRAL_NOT_FOUND_ERROR = "Cannot integrate it analytically. Please use Numerical Integration methods."
struct IntegralNotFoundError <: Exception end
function Base.showerror(io::IO, ::IntegralNotFoundError)
    return print(io, INTEGRAL_NOT_FOUND_ERROR)
end

const DERIVATIVE_NOT_FOUND_ERROR = "Derivatives greater than second order is not supported."
struct DerivativeNotFoundError <: Exception end
function Base.showerror(io::IO, ::DerivativeNotFoundError)
    return print(io, DERIVATIVE_NOT_FOUND_ERROR)
end

const INTEGRAL_INVERSE_NOT_FOUND_ERROR = "Cannot invert the integral analytically. Please use Numerical methods."
struct IntegralInverseNotFoundError <: Exception end
function Base.showerror(io::IO, ::IntegralInverseNotFoundError)
    return print(io, INTEGRAL_INVERSE_NOT_FOUND_ERROR)
end

const INTEGRAL_NOT_INVERTIBLE_ERROR = "The Interpolation is not positive everywhere so its integral is not invertible."
struct IntegralNotInvertibleError <: Exception end
function Base.showerror(io::IO, ::IntegralNotInvertibleError)
    return print(io, INTEGRAL_NOT_INVERTIBLE_ERROR)
end

const EXTRAPOLATION_NOT_IMPLEMENTED_ERROR = "The provided extrapolation option is not implemented."
struct ExtrapolationNotImplementedError <: Exception end
function Base.showerror(io::IO, ::ExtrapolationNotImplementedError)
    return print(io, EXTRAPOLATION_NOT_IMPLEMENTED_ERROR)
end

"""
    output_dim(x::AbstractInterpolation)

Return the number of dimensions `ndims(x(t))` of interpolation `x` for a scalar `t`.
"""
output_dim(x::AbstractInterpolation) = _output_dim(x.u)
_output_dim(::AbstractVector) = 0 # each value is a scalar
_output_dim(::AbstractVector{<:AbstractArray{<:Any, N}}) where {N} = N # each value is an array but values are not stacked
_output_dim(::AbstractArray{<:Any, N}) where {N} = N - 1 # each value is an array but multiple values are stacked

"""
    output_size(x::AbstractInterpolation)

Return the size `size(x(t))` of interpolation `x` for a scalar `t`.
"""
output_size(x::AbstractInterpolation) = _output_size(x.u)
_output_size(::AbstractVector{<:Number}) = ()
_output_size(u::AbstractVector) = size(first(u))
_output_size(u::AbstractArray) = Base.front(size(u))

export LinearInterpolation, QuadraticInterpolation, LagrangeInterpolation,
    AkimaInterpolation, ConstantInterpolation, SmoothedConstantInterpolation,
    QuadraticSpline, CubicSpline, BSplineInterpolation, BSplineApprox,
    CubicHermiteSpline, PCHIPInterpolation, QuinticHermiteSpline,
    SmoothArcLengthInterpolation, LinearInterpolationIntInv,
    ConstantInterpolationIntInv, ExtrapolationType
export output_dim, output_size

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
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    function RegularizationSmooth(
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
            extrapolation_right
        )
        return new{typeof(u), typeof(t), eltype(u), typeof(λ), typeof(Aitp)}(
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
            extrapolation_right
        )
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
        return new{
            typeof(u), typeof(t), typeof(m),
            typeof(p0), typeof(ub), typeof(lb),
            typeof(alg), typeof(pmin), eltype(u),
        }(
            u,
            t,
            m,
            p0,
            ub,
            lb,
            alg,
            pmin,
            extrapolate
        )
    end
end

# Define an empty function, so that it can be extended via `DataInterpolationsOptimExt`
function Curvefit()
    error("CurveFit requires loading Optim and ForwardDiff, e.g. `using Optim, ForwardDiff`")
end

export Curvefit

end # module
