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

# Build a function from the symbolic expression and evaluate it
f_ex = build_function(ex, τ, expression = Val{false})
@test f_ex(0.5) ≈ cos(0.5) * A(0.5) # true
```

### Symbolic Derivatives

```@example symbolics
D = Differential(τ)

ex1 = A(τ)

# Derivative of interpolation
ex2 = expand_derivatives(D(ex1))

# Build a function from the derivative expression and evaluate it
f_deriv = build_function(ex2, τ, expression = Val{false})
@test f_deriv(0.5) ≈ DataInterpolations.derivative(A, 0.5) # true

# Higher Order Derivatives
ex3 = expand_derivatives(D(D(A(τ))))

f_deriv2 = build_function(ex3, τ, expression = Val{false})
@test f_deriv2(0.5) ≈ DataInterpolations.derivative(A, 0.5, 2) # true
```

### Symbolic construction

The examples above evaluate a concretely-built interpolation at a symbolic *time* argument. `u` itself can also be symbolic, e.g. from `@variables u[1:n]`:

```@example symbolics
@variables u[1:5]
t2 = 0.0:1.0:4.0
B = LinearInterpolation(u, t2)
B(2.5) # a Num expression in u[1], ..., u[5]
```

Every interpolation method supports this except `SmoothArcLengthInterpolation`, whose constructor does computational geometry (circle/line segment fitting, intersection detection) that depends on the concrete shape of the data and can't be resolved symbolically.

## Using with ModelingToolkit.jl

We recommend using the
[ModelingToolkitStandardLibrary Interpolation Blocks](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/input_component/)
in order to use DataInterpolations.jl in MTK models.
