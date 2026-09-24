"""
    AbstractIntegralInverseInterpolation{T}

Developer supertype for interpolation objects that evaluate the inverse of a cumulative
integral. Subtypes retain the [`AbstractInterpolation`](@ref) call and derivative
interfaces, with `u` representing the original time values and `t` representing cumulative
integral values. The internal `itp` field stores the source interpolation used to evaluate
the inverse derivative.

Users should construct these objects with [`invert_integral`](@ref), rather than calling a
concrete constructor directly. A new implementation must provide the same
`_interpolate(A, t, iguess)` contract as [`DataInterpolations._interpolate`](@ref) and return
the source interpolation's reciprocal value from `_derivative`.

# Fields

- `u`: the original interpolation's sample locations.
- `t`: cumulative integral values used as inverse-interpolation locations.
- `itp`: the source interpolation used to evaluate the inverse derivative.

Implementations also provide the interval-search fields described for
[`AbstractInterpolation`](@ref), including `iguesser`, `kind`, and `t_props`.
"""
abstract type AbstractIntegralInverseInterpolation{T} <: AbstractInterpolation{T} end

"""
    invert_integral(A::AbstractInterpolation)::AbstractIntegralInverseInterpolation

Creates the inverted integral interpolation object from the given interpolation. Conditions:

  - The range of `A` must be strictly positive
  - `A.u` must be a number type (on which an ordering is defined)
  - This is currently only supported for `ConstantInterpolation` and `LinearInterpolation`

# Arguments

  - `A`: interpolation object satisfying the above requirements

# Returns

- An [`AbstractIntegralInverseInterpolation`](@ref) that maps an integrated value back to
  its time location.

# Examples

```julia
using DataInterpolations

A = LinearInterpolation([1.0, 2.0, 3.0], [0.0, 1.0, 2.0])
A_inverse = invert_integral(A)
A_inverse(1.0)
```
"""
invert_integral(::AbstractInterpolation) = throw(IntegralInverseNotFoundError())

_integral(::AbstractIntegralInverseInterpolation, idx, t) = throw(IntegralNotFoundError())

function _derivative(A::AbstractIntegralInverseInterpolation, t::Number, iguess)
    return inv(A.itp(_interpolate(A, t, iguess)))
end

"""
    LinearInterpolationIntInv(u, t, A)

It is the interpolation of the inverse of the integral of a `LinearInterpolation`.
Can be easily constructed with `invert_integral(A::LinearInterpolation{<:AbstractVector{<:Number}})`

# Arguments

  - `u` : Given by `A.t`
  - `t` : Given by `A.I` (the cumulative integral of `A`)
  - `A` : The `LinearInterpolation` object

# Examples

```julia
using DataInterpolations

A = LinearInterpolation([1.0, 2.0, 3.0], [0.0, 1.0, 2.0])
A_inverse = invert_integral(A)
A_inverse(1.0)
```
"""
struct LinearInterpolationIntInv{uType, tType, itpType, T, propsType} <:
    AbstractIntegralInverseInterpolation{T}
    u::uType
    t::tType
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    iguesser::Guesser{tType}
    t_props::propsType
    kind::FindFirstFunctions.StrategyKind
    itp::itpType
    function LinearInterpolationIntInv(
            u, t, A, extrapolation_left, extrapolation_right, t_props,
        )
        kind = _resolve_strategy_kind(t, t_props)
        return new{
            typeof(u), typeof(t), typeof(A), eltype(u),
            typeof(t_props),
        }(
            u, t, extrapolation_left, extrapolation_right,
            Guesser(t), t_props, kind, A
        )
    end
end

function invertible_integral(A::LinearInterpolation{<:AbstractVector{<:Number}})
    return all(A.u .> 0)
end

function get_I(A::AbstractInterpolation)
    I = isempty(A.I) ? cumulative_integral(A, true) : copy(A.I)
    pushfirst!(I, 0)
    return I
end

function invert_integral(
        A::LinearInterpolation{<:AbstractVector{<:Number}};
        extrapolation_left::ExtrapolationType.T = A.extrapolation_left,
        extrapolation_right::ExtrapolationType.T = A.extrapolation_right,
        search_properties::Union{Nothing, FindFirstFunctions.SearchProperties} = nothing
    )
    !invertible_integral(A) && throw(IntegralNotInvertibleError())
    t_I = get_I(A)
    t_props = something(search_properties, FindFirstFunctions.SearchProperties(t_I))
    return LinearInterpolationIntInv(
        A.t, t_I, A, extrapolation_left, extrapolation_right, t_props
    )
end

function _interpolate(
        A::LinearInterpolationIntInv{<:AbstractVector{<:Number}}, t::Number, iguess
    )
    idx = get_idx(A, t, iguess)
    Δt = t - A.t[idx]
    x = A.itp.u[idx]
    slope = get_parameters(A.itp, idx)
    u = A.u[idx] + 2Δt / (x + sqrt(x^2 + slope * 2Δt))
    return u
end

"""
    ConstantInterpolationIntInv(u, t, A)

It is the interpolation of the inverse of the integral of a `ConstantInterpolation`.
Can be easily constructed with `invert_integral(A::ConstantInterpolation{<:AbstractVector{<:Number}})`

# Arguments

  - `u` : Given by `A.t`
  - `t` : Given by `A.I` (the cumulative integral of `A`)
  - `A` : The `ConstantInterpolation` object

# Examples

```julia
using DataInterpolations

A = ConstantInterpolation([1.0, 2.0, 3.0], [0.0, 1.0, 2.0])
A_inverse = invert_integral(A)
A_inverse(1.0)
```
"""
struct ConstantInterpolationIntInv{uType, tType, itpType, T, propsType} <:
    AbstractIntegralInverseInterpolation{T}
    u::uType
    t::tType
    extrapolation_left::ExtrapolationType.T
    extrapolation_right::ExtrapolationType.T
    iguesser::Guesser{tType}
    t_props::propsType
    kind::FindFirstFunctions.StrategyKind
    itp::itpType
    function ConstantInterpolationIntInv(
            u, t, A, extrapolation_left, extrapolation_right, t_props,
        )
        kind = _resolve_strategy_kind(t, t_props)
        return new{
            typeof(u), typeof(t), typeof(A), eltype(u),
            typeof(t_props),
        }(
            u, t, extrapolation_left, extrapolation_right,
            Guesser(t), t_props, kind, A
        )
    end
end

function invertible_integral(A::ConstantInterpolation{<:AbstractVector{<:Number}})
    return all(A.u .> 0)
end

function invert_integral(
        A::ConstantInterpolation{<:AbstractVector{<:Number}};
        extrapolation_left::ExtrapolationType.T = A.extrapolation_left,
        extrapolation_right::ExtrapolationType.T = A.extrapolation_right,
        search_properties::Union{Nothing, FindFirstFunctions.SearchProperties} = nothing
    )
    !invertible_integral(A) && throw(IntegralNotInvertibleError())
    t_I = get_I(A)
    t_props = something(search_properties, FindFirstFunctions.SearchProperties(t_I))
    return ConstantInterpolationIntInv(
        A.t, t_I, A, extrapolation_left, extrapolation_right, t_props
    )
end

function _interpolate(
        A::ConstantInterpolationIntInv{<:AbstractVector{<:Number}}, t::Number, iguess
    )
    idx = get_idx(A, t, iguess; ub_shift = 0)
    if A.itp.dir === :left
        # :left means that value to the left is used for interpolation
        idx_ = get_idx(A, t, idx; lb = 1, ub_shift = 0)
    else
        # :right means that value to the right is used for interpolation
        idx_ = get_idx(A, t, idx; side = :first, lb = 1, ub_shift = 0)
    end
    return A.u[idx] + (t - A.t[idx]) / A.itp.u[idx_]
end
