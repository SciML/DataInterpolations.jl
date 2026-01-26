# Clarifications

Here are some numerical clarifications of the algorithm definitions and behaviors.

# 0-order B-Spline Interpolations are not the same as Constant Interpolations

```@example interpclarity
using DataInterpolations, Plots
using DataInterpolations: derivative
u = [2.0, 1.0, 5.0, 4.0, 5.0, 4.0]
t = [0.0, 2.0, 3.5, 4.0, 5.0, 6.5]
bspline = BSplineInterpolation(u, t, 0, :Uniform, :Uniform; extrapolation_left=ExtrapolationType.Extension,  extrapolation_right=ExtrapolationType.Extension)
plot(bspline)
```

 A B-spline curve is constructed using control coefficients, hence, the jump locations are at knot vectors and do not coincide with data points(which is the case for constant-interpolation).
Thus, the plot for B-Spline interpolation does not appear the same as the plot for Constant Interpolation.

# Derivative behavior of quadratic B-Spline 

```@example interpclarity
bspline = BSplineInterpolation(u, t, 2, :Uniform, :Uniform; extrapolation_left=ExtrapolationType.Extension,  extrapolation_right=ExtrapolationType.Extension)
plot(t->derivative(bspline, t))
```

The derivative becomes piecewise linear with jumps, and extrapolation can introduce sharp artifacts near boundaries, additionally derivative(bspline, t) is calculated by differentiating the spline basis and evaluating outside nominal knot spans using extension, not clamping.
Hence, sudden spikes near the ends and sharp negative excursions where basis support ends.
