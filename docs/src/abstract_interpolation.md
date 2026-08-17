# Abstract interpolation interface

`AbstractInterpolation` is the extension point for adding a new interpolation method. The
public evaluation, differentiation, and integration functions are generic over this type,
so downstream code should use those functions rather than dispatching on a concrete
interpolation implementation.

```@docs
DataInterpolations.AbstractInterpolation
```

## Required data

A concrete subtype must provide the following fields:

- `u`: data values. A vector stores scalar values or one array per sample; an array stores
  samples along its last dimension.
- `t`: ordered sample locations with `length(t)` equal to the number of samples in `u`.
- `iguesser`: reusable interval-search state used by scalar evaluation and differentiation.
- `extrapolation_left` and `extrapolation_right`: values from
  [`ExtrapolationType`](@ref) controlling out-of-range evaluation.

For the generic integral interface, the subtype also provides:

- `I`: cumulative integrals for the intervals, or an empty cache when
  `cache_parameters` is `false`. A method without an analytic integral should omit this
  field; `integral` then reports `IntegralNotFoundError`.
- `cache_parameters::Bool`: whether the generic integral method reads the cached `I` values.
- `kind` and `t_props`: interval-search state, normally constructed from
  `FindFirstFunctions.SearchProperties(t)` and `FindFirstFunctions.Auto(t, t_props)`.

Construct `iguesser` with `FindFirstFunctions.Guesser(t)`. The `t` values must remain ordered,
and `u` must store one sample per `t` entry along its last dimension.

## Developer hooks

The following methods are the implementation hooks used by the generic interface. They are
public developer API for packages that add interpolation subtypes, but they are not intended
to be called by end users.

```@docs
DataInterpolations._interpolate
DataInterpolations._derivative
DataInterpolations._integral
```

The hooks must satisfy these rules:

- `_interpolate(A, t, iguess)` evaluates one scalar location and returns one data point.
- `_derivative(A, t, iguess)` returns the first derivative with the same shape as one data
  point. At a discontinuous knot it returns the left derivative.
- `_integral(A, idx, t1, t2)` integrates the requested subinterval of
  `[A.t[idx], A.t[idx + 1]]` and preserves the data-point shape.

The generic forms `A(t)`, `derivative(A, t, order = 1)`, and `integral(A, t)` compose these
hooks with interval search and extrapolation. A method should therefore be tested through
those generic functions, including scalar, vectorized, and in-place evaluation where its
data shape supports them.

## Minimal implementation

The following sketch shows the required dispatch. A production implementation should also
validate its data and provide the fields needed by any optional features it supports.

```julia
using FindFirstFunctions

struct MyInterpolation <: DataInterpolations.AbstractInterpolation{Float64}
    u::Vector{Float64}
    t::Vector{Float64}
    I::Vector{Float64}
    iguesser
    extrapolation_left::DataInterpolations.ExtrapolationType.T
    extrapolation_right::DataInterpolations.ExtrapolationType.T
    kind::FindFirstFunctions.StrategyKind
    t_props::FindFirstFunctions.SearchProperties
    cache_parameters::Bool
end

DataInterpolations._interpolate(A::MyInterpolation, t::Number, iguess) = 2t
DataInterpolations._derivative(A::MyInterpolation, t::Number, iguess) = 2.0
DataInterpolations._integral(
    A::MyInterpolation, idx::Integer, t1::Number, t2::Number
) = 2 * (t2 - t1)
```
