# Interface for using the Interpolations object

We will again use the same data as the previous tutorial to demonstrate how to use the Interpolations object for computing interpolated values at any time point, its derivatives and integrals.

```@example interface
using DataInterpolations

# Dependent variable
u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]

# Independent variable
t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]
```

## Interpolated values

All Interpolation methods return an object from which we can compute the value of the dependent variable at any time point.

We will use `CubicSpline` method for demonstration but the API is same for all the methods.

```@example interface
A = CubicSpline(u, t)

# For interpolation do, A(t)
A(100.0)
```

!!! note
    
    The values computed beyond the range of the time points provided during interpolation will not be reliable as these methods only perform well within the range and the first/last piece polynomial fit is extrapolated on either sides which might not reflect the true nature of the data.

## Derivatives

Derivatives of the interpolated curves can also be computed at any point for all the methods.

We will continue with the above example, but the API is same for all the methods.

```@example interface
# derivative(A, t)
DataInterpolations.derivative(A, 1.0)
```

## Integrals

Integrals of the interpolated curves can also be computed easily.

Currently, this is implemented only for a few methods - `ConstantInterpolation`, `LinearInterpolation`, `QuadraticInterpolation`, `QuadraticSpline` and `CubicSpline`.

To compute the integrals from the start of time points provided during interpolation to any point, we can do:

```@example interface
# integral(A, t)
DataInterpolations.integral(A, 5.0)
```

If we want to compute integrals between two points, we can do:

```@example interface
# integral(A, t1, t2)
DataInterpolations.integral(A, 1.0, 5.0)
```

!!! note
    
    If the times provided in the integral goes beyond the range of the time points provided during interpolation, it uses extrapolation methods to compute the values and hence the integral can be misrepsentative and might not reflect the true nature of the data.
