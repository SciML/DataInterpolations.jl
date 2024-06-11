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
affiliations:
 - name: JuliaHub
   index: 1
 - name: Massachusetts Institute of Technology
   index: 2
 - name: Pumas-AI
   index: 3
date: 6 June 2024
bibliography: paper.bib
---

# Summary

Interpolations are used to estimate values between known data points using an approximate continuous function.DataInterpolations.jl is a Julia [@Bezanson2017] package containing 1D implementations of some of the most commonly used interpolation functions. These include Constant Interpolation, Linear Interpolation, Quadratic Interpolation, Lagrange Interpolation [@lagrange], Quadratic Splines, Cubic Splines [@Schoenberg1988], Akima Splines [@10.1145/321607.321609], B-Splines [@Curry1988] [@DEBOOR197250] and Regression based B-Splines. Along with these, the package also has methods to fit parameterized curves with the data points and Tikhonov regularization [@Tikhonov1943OnTS] [@amt-14-7909-2021] for obtaining smooth curves. The package also provides functionality to compute integrals and derivatives upto second order for those interpolations methods.

# Statement of need

Interpolations are a very important component of many modeling workflows. In many models, inputs which are sampled or measured need to be represented as a continuous function or a smooth curve for simulation. In many scientific machine learning workflows, we need interpolations of data to learn continuous models. There already have been a few interpolation packages in Julia like Interpolations.jl but it has a limitation of assuming uniformly spaced data which is not usually the case with data collected from real world. DataInterpolations.jl provides fast interpolation methods for arbitrary spaced 1D data with a consistent and simple interface. It is also automatic differentiation friendly. It can also be used symbolically with Symbolics.jl [@gowda2021high] and plugged into models defined using ModelingToolkit.jl [@ma2021modelingtoolkit].

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
