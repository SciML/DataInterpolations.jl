# DataInterpolations v5 Release Notes

## Breaking changes

  - `AbstractInterpolation` is not a subtype of `AbstractVector` anymore. This was needed for previous versions of ModelingToolkit.jl to represent splines as vectors.

  - Indexing overloads for `AbstractInterpolation` and the type parameter associated with it are removed. For example - `A` is an interpolation object:
    
      + Doing `A[i]` will error. Use `A.u[i]`.
      + `size(A)` will error. Use `size(A.u)` or `size(A.t)`.
  - Removed deprecated bindings for `ZeroSpline` which is the same as `ConstantInterpolation`.
