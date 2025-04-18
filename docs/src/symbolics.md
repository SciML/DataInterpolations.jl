# Using DataInterpolations.jl with Symbolics.jl and ModelingToolkit.jl

All interpolation methods can be integrated with [Symbolics.jl](https://symbolics.juliasymbolics.org/stable/) and [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) seamlessly.

## Using with Symbolics.jl

### Expressions

```@example symbolics
using DataInterpolations, Symbolics
using Test

u = [0.0, 1.5, 0.0]
t = [0.0, 0.5, 1.0]
A = LinearInterpolation(u, t)

@variables τ

# Simple Expression
ex = cos(τ) * A(τ)
@test substitute(ex, Dict(τ => 0.5)) == cos(0.5) * A(0.5) # true
```

### Symbolic Derivatives

```@example symbolics
D = Differential(τ)

ex1 = A(τ)

# Derivative of interpolation
ex2 = expand_derivatives(D(ex1))

@test substitute(ex2, Dict(τ => 0.5)) == DataInterpolations.derivative(A, 0.5) # true

# Higher Order Derivatives
ex3 = expand_derivatives(D(D(A(τ))))

@test substitute(ex3, Dict(τ => 0.5)) == DataInterpolations.derivative(A, 0.5, 2) # true
```

## Using with ModelingToolkit.jl

We recommend using the
[ModelingToolkitStandardLibrary Interpolation Blocks](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/input_component/)
in order to use DataInterpolations.jl in MTK models.
