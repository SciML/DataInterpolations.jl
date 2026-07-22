module DataInterpolations

### Interface Functionality

abstract type AbstractInterpolation{T} end

using LinearAlgebra: LinearAlgebra, Tridiagonal, dot, norm, normalize!
using RecipesBase: RecipesBase, @recipe, @series
using PrettyTables: PrettyTables, pretty_table
using ForwardDiff: ForwardDiff
using EnumX: EnumX, @enumx
using StaticArrays: SVector
import FindFirstFunctions
import FindFirstFunctions: Guesser

"""
    ExtrapolationType

Enumeration of extrapolation modes used by interpolation constructors through
the `extrapolation`, `extrapolation_left`, and `extrapolation_right` keyword
arguments.

## Values

  - `ExtrapolationType.None`: throw an error outside the interpolation interval.
  - `ExtrapolationType.Constant`: use the nearest endpoint value.
  - `ExtrapolationType.Linear`: extend the interpolation linearly from the nearest interval.
  - `ExtrapolationType.Extension`: evaluate the interpolation formula outside the data range.
  - `ExtrapolationType.Periodic`: wrap the query point periodically onto the data range.
  - `ExtrapolationType.Reflective`: reflect the query point back onto the data range.

## Examples

```julia
using DataInterpolations

t = [0.0, 1.0]
u = [1.0, 3.0]
A = LinearInterpolation(u, t; extrapolation = ExtrapolationType.Extension)

A(1.5)
```
"""
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

"""
    Curvefit(
            u, t, model, p0, alg, box = false, lb = nothing, ub = nothing;
            extrapolate = false
        ) -> CurvefitCache

Fit `model(t, p)` to data values `u` at sample locations `t` with Optim.jl and return the
resulting callable interpolation. The returned interpolation evaluates `model` with the
optimized parameters available as its `pmin` field.

# Arguments

- `u`: observed scalar data values.
- `t`: sample locations corresponding to `u`.
- `model`: callable with signature `model(t, p)`, where `p` is a parameter vector.
- `p0`: initial parameter vector supplied to Optim.jl.
- `alg`: Optim.jl optimization algorithm, such as `LBFGS()`.
- `box::Bool = false`: use boxed optimization when `true`; both bounds must then be given.
- `lb = nothing`: lower parameter bounds used when `box` is `true`.
- `ub = nothing`: upper parameter bounds used when `box` is `true`.

# Keywords

- `extrapolate::Bool = false`: permit evaluation outside the range of `t`. When `false`,
  out-of-range evaluation throws `ExtrapolationError`.

# Returns

- `CurvefitCache`: a callable fitted interpolation. Its `pmin` field contains the optimized
  parameter vector.

# Examples

```julia
using DataInterpolations, Optim

model(t, p) = p[1] .* t .+ p[2]
A = Curvefit([1.0, 3.0, 5.0], [0.0, 1.0, 2.0], model, [0.0, 0.0], LBFGS())

A(1.5)
```
"""
function Curvefit()
    error("CurveFit requires loading Optim and ForwardDiff, e.g. `using Optim, ForwardDiff`")
end

export Curvefit

end # module
