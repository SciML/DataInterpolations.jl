# Interface for using the Interpolations object

We will again use the same data as the previous tutorial to demonstrate how to use the Interpolations object for computing interpolated values at any time point,
as well as derivatives and integrals.

```@example interface
using DataInterpolations

# Dependent variable
u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]

# Independent variable
t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
```

## Interpolated values

All interpolation methods return an object from which we can compute the value of the dependent variable at any time point.

We will use the `CubicSpline` method for demonstration, but the API is the same for all the methods. We can also pass the `extrapolation = ExtrapolationType.Extension` keyword if we want to allow the interpolation to go beyond the range of the timepoints in the positive `t` direction. The default value is `extrapolation = ExtrapolationType.None`. For more information on extrapolation see [Extrapolation methods](extrapolation_methods.md).

```@example interface
A1 = CubicSpline(u, t)

# For interpolation do, A(t)
A1(100.0)

A2 = CubicSpline(u, t; extrapolation = ExtrapolationType.Extension)

# Extrapolation
A2(300.0)
```

!!! note

    The values computed beyond the range of the time points provided during interpolation will not be reliable, as these methods only perform well within the range and the first/last piece polynomial fit is extrapolated on either side which might not reflect the true nature of the data.

The keyword `cache_parameters = true` can be passed to precalculate parameters at initialization, making evaluations cheaper to compute. This is not compatible with modifying `u` and `t`. The default `cache_parameters = false` does however not prevent allocation in every interpolation constructor call.

### In-place Evaluation (Allocation-free)

When performance is critical, such as in ODE solvers or tight loops, you can use the in-place variant to avoid memory allocations. This is particularly useful when you need to interpolate at many points repeatedly.

To use in-place interpolation, pass a pre-allocated output array as the first argument:

```@example interface
# Pre-allocate output array
t_eval = [50.0, 100.0, 150.0, 200.0]
u_out = similar(u, length(t_eval))

# In-place interpolation: interp(output, t_values)
A1(u_out, t_eval)

u_out
```

The in-place form `interp(out, t)` writes the interpolated values directly into `out`, avoiding allocation of a new array. The output array must have the same length as the input time vector.

For vector-of-arrays data (where `u` is a `Vector{Vector}` or `Vector{Matrix}`), the output should be a pre-allocated vector of arrays with the same structure:

```@example interface
# Vector of vectors example
u_vov = [[14.7, 7.35], [11.51, 5.76], [10.41, 5.21], [14.95, 7.48], [12.24, 6.12], [11.22, 5.61]]
A_vov = CubicSpline(u_vov, t)

# Pre-allocate output as Vector{Vector}
t_eval_3 = [50.0, 100.0, 150.0]
out_vov = [zeros(2) for _ in 1:length(t_eval_3)]
A_vov(out_vov, t_eval_3)

out_vov
```

For multi-dimensional data (where `u` is a matrix or higher-dimensional array), the output array should have the same leading dimensions as `u`, with the last dimension matching the length of `t`:

```@example interface
# Matrix example (stacked form)
u_matrix = [14.7 11.51 10.41 14.95 12.24 11.22;
            7.35 5.76 5.21 7.48 6.12 5.61]
A_matrix = CubicSpline(u_matrix, t)

# Pre-allocate for 3 evaluation points
out_matrix = zeros(2, 3)
A_matrix(out_matrix, t_eval_3)

out_matrix
```

!!! tip "Performance tip"

    Using in-place interpolation can significantly reduce allocations in performance-critical code. For example, in an ODE derivative function, switching from `interp(t_vec)` to `interp(out, t_vec)` can eliminate allocations entirely within the hot loop.

## Derivatives

Derivatives of the interpolated curves can also be computed at any point for all the methods. Derivatives up to second order are supported where first order derivative is computed analytically and second order using `ForwardDiff.jl`. Order is passed as the third argument. It is 1 by default.

We will continue with the above example, but the API is the same for all the methods. If the interpolation is defined with `extrapolate=true`, derivatives can also be extrapolated.

```@example interface
# derivative(A, t)
DataInterpolations.derivative(A1, 1.0, 1)
DataInterpolations.derivative(A1, 1.0, 2)

# Extrapolation
DataInterpolations.derivative(A2, 300.0)
```

## Integrals

Integrals of the interpolated curves can also be computed easily.

!!! note
    
    Integrals for `LagrangeInterpolation`, `BSplineInterpolation`, `BSplineApprox`, `Curvefit` will error as there are no simple analytical solutions available. Please use numerical methods instead, such as [Integrals.jl](https://docs.sciml.ai/Integrals/stable/).

To compute the integrals from the start of time points provided during interpolation to any point, we can do:

```@example interface
# integral(A, t)
DataInterpolations.integral(A1, 5.0)
```

If we want to compute integrals between two points, we can do:

```@example interface
# integral(A, t1, t2)
DataInterpolations.integral(A1, 1.0, 5.0)
```

Again, if the interpolation is defined with `extrapolate=true`, the integral can be computed beyond the range of the timepoints.

```@example interface
# integral(A, t1, t2)
DataInterpolations.integral(A2, 200.0, 300.0)
```

!!! note

    If the times provided in the integral go beyond the range of the time points provided during interpolation, it uses extrapolation methods to compute the values, and hence the integral can be misrepresentative and might not reflect the true nature of the data.
