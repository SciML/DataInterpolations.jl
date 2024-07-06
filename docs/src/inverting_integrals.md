# Inverting integrals

Solving implicit integral problems of the form:

```math
\begin{equation}
    \text{find $t$ such that } \int_{t_1}^t f(\tau)\text{d}\tau = V \ge 0
\end{equation}
```

is supported for interpolations $f$ that are strictly positive and of one of these types:

  - `ConstantInterpolation`
  - `LinearInterpolation`

This is done by creating an 'integral inverse' interpolation object which can efficiently compute $t$ for a given value of $V$, see the example below.

```@example inverting_integrals
using Random #hide
Random.seed!(1234) # hide
using DataInterpolations
using Plots

# Create LinearInterpolation object from the
u = sqrt.(1:25) + (2.0 * rand(25) .- 1.0) / 3
t = cumsum(rand(25))
A = LinearInterpolation(u, t)

# Create LinearInterpolationIntInv object
# from the LinearInterpolation object
A_intinv = DataInterpolations.invert_integral(A)

# Get the t values up to and including the
# solution to the integral problem
V = 25.0
t_ = A_intinv(V)
ts = t[t .<= t_]
push!(ts, t_)

# Plot results
plot(A; label = "Linear Interpolation")
plot!(ts, A.(ts), fillrange = 0.0, fillalpha = 0.75,
    fc = :blues, lw = 0, label = "Area of $V")
```

## Docstrings

```@docs
DataInterpolations.invert_integral
ConstantInterpolationIntInv
LinearInterpolationIntInv
```
