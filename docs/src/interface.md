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

The keyword `cache_parameters = true` can be passed to precalculate parameters at initialization, making evalations cheaper to compute. This is not compatible with modifying `u` and `t`. The default `cache_parameters = false` does however not prevent allocation in every interpolation constructor call.

## Derivatives

Derivatives of the interpolated curves can also be computed at any point for all the methods. Derivatives upto second order is supported where first order derivative is computed analytically and second order using `ForwardDiff.jl`. Order is passed as the third argument. It is 1 by default.

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
    
    If the times provided in the integral go beyond the range of the time points provided during interpolation, it uses extrapolation methods to compute the values, and hence the integral can be misrepsentative and might not reflect the true nature of the data.
