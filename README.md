# DataInterpolations.jl

DataInterpolations.jl is a library for performing interpolations of one-dimensional data. By 
"data interpolations" we mean techniques for interpolating possibly noisy data, and thus
some methods are mixtures of regressions with interpolations (i.e. do not hit the data
points exactly, smoothing out the lines). This library can be used to fill in intermediate
data points in applications like timeseries data.

## API

All interpolation objects act as functions. Thus for example, using an interpolation looks like:

```julia
u = rand(5)
t = 0:4
interp = LinearInterpolation(u,t)
interp(3.5) # Gives the linear interpolation value at t=3.5
```

The indexing retreives the underlying values:

```julia
interp[4] # Gives the 4th value of u
```

## Available Interpolations

In all cases, `u` an `AbstractVector` of values and `t` is an `AbstractVector` of timepoints
corresponding to `(t,u)` pairs.

- `ZeroSpline(u,t)` - A piecewise constant interpolation.
- `LinearInterpolation(u,t)` - A linear interpolation.
- `QuadraticInterpolation(u,t)` - A quadratic interpolation.
- `LagrangeInterpolation(u,t,n)` - A Lagrange interpolation of order `n`.
- `QuadraticSpline(u,t)` - A quadratic spline interpolation.
- `CubicSpline(u,t)` - A cubic spline interpolation.
- `BSpline(u,t,d,pVec,knotVec)` - A regression B-spline.
- `Loess(u,t,d,α)` - A local least square regression. `d` is the number of local intervals to use and `α` is a smoothing parameter.
- `GPInterpolation(u,t,m,k,n=-2.0)` - A Gaussian Process interpolation. Stochastic: each trajectory is different. `m` is ... . `k` is ... . `n` is ... and defaults to `-2.0`.
- `Curvefit(u,t,m,p)` - An interpolation which is done by fitting a user-given functional form `m(t,p)` where `p` is the vector of parameters. The user's input `p` is a an initial value for a least-square fitting and the optimal `p`s are used in the interpolation.
- `SigmoidFit(u,t,p=zeros(2))` - A `CurveFit` to the Sigmoid function: `m(t,p) = p[1]/(1+exp(t-p[2]))`. Initial parameters default to zeros.
- `HillFit(u,t,p=zeros(2),n=nothing)` - A `CurveFit` to the Hill function: `m(t,p) = p[1] * inv(1.0 + p[2]/(t^n))`. If `n` is given then it's assumed constant, otherwise it is a parameter which is fit. Initial parameters default to zero.
- `WeibullFit(u,t,p=zeros(3))` - A `CurveFit` to the Weibull exponential function: `m(t,p) = max(0,p[1] * (1 - exp(-1 * (t/p[2])^p[3])))`. Initial parameters default to zero.
