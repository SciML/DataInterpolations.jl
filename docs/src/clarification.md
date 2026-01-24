# Clarifications
We will clarify regarding the issues raised through this doc

# 0-order B-Spline Interpolations are not the same as Constant Interpolations
```@example tutorial
julia> using DataInterpolations, Plots
julia> using DataInterpolations: derivative
julia> u = [2.0, 1.0, 5.0, 4.0, 5.0, 4.0]
julia> t = [0.0, 2.0, 3.5, 4.0, 5.0, 6.5]
julia> bspline = BSplineInterpolation(u, t, 0, :Uniform, :Uniform; extrapolation_left=ExtrapolationType.Extension,  extrapolation_right=ExtrapolationType.Extension)
julia> plot(bspline)
```
 A B-spline curve is constructed using control coefficients, hence, the jump locations are at knot vectors and do not coincide with data points(which is the case for constant-interpolation).
Thus, the plot for B-Spline interpolation does not appear the same as the plot for Constant Interpolation.

# Derivative behavior of quadratic B-Spline 
```@example tutorial
bspline = BSplineInterpolation(u, t, 2, :Uniform, :Uniform; extrapolation_left=ExtrapolationType.Extension,  extrapolation_right=ExtrapolationType.Extension)
plot(t->derivative(bspline, t))
```
The derivative becomes piecewise linear with jumps, and extrapolation can introduce sharp artifacts near boundaries, additionally derivative(bspline, t) is calculated by differentiating the spline basis and evaluating outside nominal knot spans using extension, not clamping.
Hence, sudden spikes near the ends and sharp negative excursions where basis support ends.
