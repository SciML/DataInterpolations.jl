---
title: 'DataInterpolations.jl: Fast Interpolations of 1D data'
tags:
  - julia
  - interpolations
authors:
  - name: Sathvik Bhagavan
    orcid: 0000-0003-0785-3586
    corresponding: true
    affiliation: 1
  - name: Christopher Rackauckas
    orcid: 0000-0001-5850-0663
    affiliation: "1, 2, 3"
  - name: Shubham Maddhashiya
    affiliation: 3
  - name: Bart de Koning
    orcid: 0009-0005-6134-6608
    affiliation: 4
affiliations:
 - name: JuliaHub
   index: 1
 - name: Massachusetts Institute of Technology
   index: 2
 - name: Pumas-AI
   index: 3
 - name: Deltares
   index: 4
date: 6 June 2024
bibliography: paper.bib
---

# Summary

Interpolations are used to estimate values between known data points using an approximate continuous function.DataInterpolations.jl is a Julia [@Bezanson2017] package containing 1D implementations of some of the most commonly used interpolation functions. These include Constant Interpolation, Linear Interpolation, Quadratic Interpolation, Lagrange Interpolation [@lagrange], Quadratic Splines, Cubic Splines [@Schoenberg1988], Akima Splines [@10.1145/321607.321609], Cubic Hermite Splines, Piecewise Cubic Hermite Interpolating Polynomial (PCHIP), Quintic Hermite Splines, B-Splines [@Curry1988] [@DEBOOR197250] and Regression based B-Splines. Along with these, the package also has methods to fit parameterized curves with the data points and Tikhonov regularization [@Tikhonov1943OnTS] [@amt-14-7909-2021] for obtaining smooth curves. The package also provides functionality to compute integrals and derivatives upto second order for those interpolations methods. It is also automatic differentiation friendly. It can also be used symbolically with Symbolics.jl [@gowda2021high] and plugged into models defined using ModelingToolkit.jl [@ma2021modelingtoolkit].

# Statement of need

Interpolations are a very important component of many modeling workflows. Often, sampled or measured inputs need to be transformed into continuous functions or smooth curves for simulation purposes. In many scientific machine learning workflows, interpolating data is essential to learn continuous models. DataInterpolations.jl can be used for facilitating these types of workflows. Several interpolation packages already exist in Julia, such as [Interpolations.jl](https://juliamath.github.io/Interpolations.jl/stable/), which primarily specializes in B-Splines and uniformly spaced data with some support for irregularly spaced data. In contrast, DataInterpolations.jl does not assume any specific structure in the data, offering greater flexibility for diverse datasets. [Interpolations.jl](https://juliamath.github.io/Interpolations.jl/stable/) also doesn't offer methods like Quadratic Interpolation, Lagrange Interpolation, Hermite Splines etc. [BasicInterpolators.jl](https://github.com/markmbaum/BasicInterpolators.jl) is more similar to DataInterpolations.jl, although it doesn't offer methods like B-Splines. Rest of the interpolation packages focus on particular methods like [BSplineKit.jl](https://github.com/jipolanco/BSplineKit.jl) for B-Splines, [FastChebInterp.jl](https://github.com/JuliaMath/FastChebInterp.jl) for Chebyshev interpolation, [PCHIPInterpolation](https://github.com/gerlero/PCHIPInterpolation.jl) for PCHIP interpolation etc. In summary, DataInterpolations.jl is more generic from other packages and offers many fast interpolation methods for arbitrarily spaced 1D data, all within a consistent and simple interface.

# Example

The following tutorials in the documentation [1](https://docs.sciml.ai/DataInterpolations/stable/methods/) provides how to define each of the interpolation methods and compute the value at any point. [2](https://docs.sciml.ai/DataInterpolations/stable/interface/) provides explanation for using the interface and interpolated objects for evaluating at any point, computing the derivative at any point and computing the integral between any two points.

A simple demonstration here:

```julia
using DataInterpolations

# Dependent variable
u = [14.7, 11.51, 10.41, 14.95, 12.24, 11.22]

# Independent variable
t = [0.0, 62.25, 109.66, 162.66, 205.8, 252.3]

A1 = CubicSpline(u, t)

# For interpolation do, A(t)
A1(100.0)

# derivative
## first order
DataInterpolations.derivative(A1, 1.0, 1)

## second order
DataInterpolations.derivative(A1, 1.0, 2)

# integral
DataInterpolations.integral(A1, 1.0, 5.0)
```

# References
