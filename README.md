# DataInterpolations.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/DataInterpolations/stable/)

[![codecov](https://codecov.io/gh/SciML/DataInterpolations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/DataInterpolations.jl)
[![CI](https://github.com/SciML/DataInterpolations.jl/actions/workflows/Tests.yml/badge.svg?branch=master)](https://github.com/SciML/DataInterpolations.jl/actions/workflows/Tests.yml)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.06917/status.svg)](https://doi.org/10.21105/joss.06917)

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
interp = LinearInterpolation(u, t)
interp(3.5) # Gives the linear interpolation value at t=3.5
```

We can efficiently interpolate onto a vector of new `t` values:

```julia
t′ = 0.5:1.0:3.5
interp(t′)
```

In-place interpolation also works:

```julia
u′ = similar(u, length(t′))
interp(u′, t′)
```

## Available Interpolations

In all cases, `u` an `AbstractVector` of values and `t` is an `AbstractVector` of timepoints
corresponding to `(u,t)` pairs.

  - `ConstantInterpolation(u,t)` - A piecewise constant interpolation.

  - `LinearInterpolation(u,t)` - A linear interpolation.
  - `QuadraticInterpolation(u,t)` - A quadratic interpolation.
  - `LagrangeInterpolation(u,t,n)` - A Lagrange interpolation of order `n`.
  - `QuadraticSpline(u,t)` - A quadratic spline interpolation.
  - `CubicSpline(u,t)` - A cubic spline interpolation.
  - `AkimaInterpolation(u, t)` - Akima spline interpolation provides a smoothing effect and is computationally efficient.
  - `BSplineInterpolation(u,t,d,pVec,knotVec)` - An interpolation B-spline. This is a B-spline which hits each of the data points. The argument choices are:
    
      + `d` - degree of B-spline
      + `pVec` - Symbol to Parameters Vector, `pVec = :Uniform` for uniform spaced parameters and `pVec = :ArcLen` for parameters generated by chord length method.
      + `knotVec` - Symbol to Knot Vector, `knotVec = :Uniform` for uniform knot vector, `knotVec = :Average` for average spaced knot vector.
  - `BSplineApprox(u,t,d,h,pVec,knotVec)` - A regression B-spline which smooths the fitting curve. The argument choices are the same as the `BSplineInterpolation`, with the additional parameter `h<length(t)` which is the number of control points to use, with smaller `h` indicating more smoothing.
  - `CubicHermiteSpline(du, u, t)` - A third order Hermite interpolation, which matches the values and first (`du`) order derivatives in the data points exactly.
  - `PCHIPInterpolation(u, t)` - a type of `CubicHermiteSpline` where the derivative values `du` are derived from the input data in such a way that the interpolation never overshoots the data.
  - `QuinticHermiteSpline(ddu, du, u, t)` - A fifth order Hermite interpolation, which matches the values and first (`du`) and second (`ddu`) order derivatives in the data points exactly.

## Extension Methods

The follow methods require extra dependencies and will be loaded as package extensions.

  - `Curvefit(u,t,m,p,alg)` - An interpolation which is done by fitting a user-given functional form `m(t,p)` where `p` is the vector of parameters. The user's input `p` is a an initial value for a least-square fitting, `alg` is the algorithm choice to use for optimize the cost function (sum of squared deviations) via `Optim.jl` and optimal `p`s are used in the interpolation. Requires `using Optim`.
  - `RegularizationSmooth(u,t,d;λ,alg)` - A regularization algorithm (ridge regression) which is done by minimizing an objective function (l2 loss + derivatives of order `d`) integrated in the time span. It is a global method and creates a smooth curve.
    Requires `using RegularizationTools`.

## Plotting

DataInterpolations.jl is tied into the Plots.jl ecosystem, by way of RecipesBase. Any interpolation can be plotted using the `plot` command (or any other), since they have type recipes associated with them. For convenience, and to allow keyword arguments to propagate properly, DataInterpolations.jl also defines several series types, corresponding to different interpolations.

The series types defined are:

  - `:linear_interp`
  - `:quadratic_interp`
  - `:lagrange_interp`
  - `:quadratic_spline`
  - `:cubic_spline`
  - `:akima_interp`
  - `:bspline_interp`
  - `:bspline_approx`
  - `:cubic_hermite_spline`
  - `:pchip_interp`
  - `:quintic_hermite_spline`

By and large, these accept the same keywords as their function counterparts.

Some keywords differ from regular plots. `label_interp` is used to label the interpolation line plot, while `label_data` labels the data points. By default, both are plotted in the same color.

## Citing

If you use this software in your work, please cite:

```bib
@article{Bhagavan2024,
  doi = {10.21105/joss.06917},
  url = {https://doi.org/10.21105/joss.06917},
  year = {2024},
  publisher = {The Open Journal},
  volume = {9},
  number = {101},
  pages = {6917},
  author = {Sathvik Bhagavan and Bart de Koning and Shubham Maddhashiya and Christopher Rackauckas},
  title = {DataInterpolations.jl: Fast Interpolations of 1D data},
  journal = {Journal of Open Source Software}
}
```
