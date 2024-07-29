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

Most common use case with [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) is to plug in interpolation objects as input functions. This can be done using `TimeVaryingFunction` component of [ModelingToolkitStandardLibrary.jl](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/).

```@example mtk
using DataInterpolations
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq

us = [0.0, 1.5, 0.0]
times = [0.0, 0.5, 1.0]
A = LinearInterpolation(us, times)

@named src = TimeVaryingFunction(A)
vars = @variables x(t) out(t)
eqs = [out ~ src.output.u, D(x) ~ 1 + out]
@named sys = ODESystem(eqs, t, vars, []; systems = [src])

sys = structural_simplify(sys)
prob = ODEProblem(sys, [x => 0.0], (times[1], times[end]))
sol = solve(prob)
```
