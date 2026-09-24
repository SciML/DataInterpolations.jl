# DataInterpolations v9 Release Notes

## Breaking changes

  - The deprecated `RegularizationTools` extension and the `RegularizationSmooth` interpolation type have been removed. `RegularizationTools` was deprecated and capped `Optim` to `≤ 1`; removing it restores support for `Optim` v2.

  - The `assume_linear_t` constructor keyword and the `DataInterpolations.looks_linear` utility have been removed. Knot-vector structure is now probed once at construction through `FindFirstFunctions.SearchProperties(t)` and cached on every interpolation as `A.t_props`; uniformly-spaced knots are detected exactly and automatically. To override the probe, pass the new `search_properties` keyword accepted by every constructor, e.g. `LinearInterpolation(u, t; search_properties = FindFirstFunctions.SearchProperties(t; is_uniform = true))`.

## New features

  - Every interpolation constructor accepts a `search_properties::Union{Nothing, FindFirstFunctions.SearchProperties}` keyword. The default `nothing` probes `t` at construction; passing a pre-built `SearchProperties` skips the probe (useful when constructing many interpolations over the same knot vector).

  - Knot search is dispatched through `FindFirstFunctions.Auto(t)` resolved at construction: uniformly-spaced knots (any `AbstractRange`, or vectors detected as exactly uniform) use a closed-form O(1) lookup; short non-uniform knot vectors use a linear scan; everything else keeps the previous bracketed gallop.

  - `LinearInterpolation` with uniformly-spaced knots and floating-point values takes a statically-dispatched fast path (closed-form index + lerp, verified against the live knots) — 5-10x faster per query on uniform grids.

  - `QuadraticSpline` construction is now O(n) instead of O(n^2) (running locator in `quadratic_spline_params`), e.g. ~870x faster at 100k points.

# DataInterpolations v5 Release Notes

## Breaking changes

  - `AbstractInterpolation` is not a subtype of `AbstractVector` anymore. This was needed for previous versions of ModelingToolkit.jl to represent splines as vectors.

  - Indexing overloads for `AbstractInterpolation` and the type parameter associated with it are removed. For example - `A` is an interpolation object:
    
      + Doing `A[i]` will error. Use `A.u[i]`.
      + `size(A)` will error. Use `size(A.u)` or `size(A.t)`.
  - Removed deprecated bindings for `ZeroSpline` which is the same as `ConstantInterpolation`.

# DataInterpolations v6 Release Notes

## Breaking changes

  - https://github.com/SciML/DataInterpolations.jl/pull/274 introduced caching of parameters for interpolations (released in v5.3) and also introduced a field `safetycopy` which was a boolean flag to create a copy of the data as the parameters would be invalid if data is mutated. This was removed in https://github.com/SciML/DataInterpolations.jl/pull/315 to introduce `cache_parameters` which made it explicit if a user wants to opt in for parameter caching or not.
