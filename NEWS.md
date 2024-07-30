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
